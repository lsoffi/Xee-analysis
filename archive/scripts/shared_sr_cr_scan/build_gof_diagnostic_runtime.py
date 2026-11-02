"""Build an isolated diagnostic-only Combine runtime using installed offline inputs.

Only GoodnessOfFit.cc is recompiled. Other object files are copied and hashed.
The production source, build directory, environment and executables are read-only.
"""
from pathlib import Path
import difflib
import hashlib
import json
import shlex
import shutil
import subprocess
import time

PROJECT = Path(__file__).resolve().parents[2]
ORIGINAL = Path('/Users/imb/Documents/Codex/2026-10-27/referenced-chatgpt-conversation-this-is-an/work/CombinedLimit')
BUILD = ORIGINAL.parent / 'build_combine'
ENV = ORIGINAL / '.pixi/envs/default'
DEST = PROJECT / 'work/shared_gof_diagnostic_runtime/v2'

CAPTURE = r'''
  // Optional diagnostics only: record the NOMINAL state before saturated fitting.
  // Values live in map nodes, whose addresses remain stable for TTree branches.
  // No additional minimization, Hesse, parameter reset or likelihood term is used.
  if (runtimedef::get("GOF_TRACK_NOMINAL")) {
    static std::map<std::string, double> nominalDiagnostics;
    auto record = [&](const std::string &name, double value) {
      const std::string branch = "gof_nominal_" + name;
      auto inserted = nominalDiagnostics.emplace(branch, value);
      if (inserted.second) {
        Combine::addBranch(branch.c_str(), &inserted.first->second,
                           (branch + "/D").c_str());
      } else {
        inserted.first->second = value;
      }
    };
    record("minimize_ok", nominalMinimizeOk ? 1. : 0.);
    record("nll", nll_nominal);
    for (RooAbsArg *arg : w->allVars()) {
      const std::string name = arg->GetName();
      const bool isNuisance = mc_s->GetNuisanceParameters() &&
                             mc_s->GetNuisanceParameters()->find(name.c_str());
      if (name == "r" || name.rfind("local_", 0) == 0 || isNuisance ||
          name.rfind("CMS_", 0) == 0 || name.rfind("Ymean_", 0) == 0) {
        if (auto *value = dynamic_cast<RooRealVar *>(arg)) record(name, value->getVal());
      }
    }
    for (RooAbsArg *arg : w->allCats()) {
      const std::string name = arg->GetName();
      if (name.rfind("pdfindex_", 0) == 0) {
        if (auto *category = dynamic_cast<RooCategory *>(arg)) record(name, category->getIndex());
      }
    }
    // save() serializes the existing minimizer state; it does not fit or run Hesse.
    // freezeDisassociatedParams performs a final fit of the selected family in
    // this runtime. Also record an explicit numerical-state consistency check.
    std::unique_ptr<RooFitResult> nominalResult(minimn.save());
    record("fitresult_available", nominalResult ? 1. : 0.);
    record("status", nominalResult ? nominalResult->status() : -999.);
    record("edm", nominalResult ? nominalResult->edm() : std::numeric_limits<double>::quiet_NaN());
    record("covQual", nominalResult ? nominalResult->covQual() : -999.);
    double maxDifference = 0.;
    if (nominalResult) {
      for (RooAbsArg *arg : nominalResult->floatParsFinal()) {
        auto *value = dynamic_cast<RooRealVar *>(arg);
        auto *current = w->var(arg->GetName());
        if (value && current) maxDifference = std::max(maxDifference, std::abs(value->getVal() - current->getVal()));
      }
    } else {
      maxDifference = std::numeric_limits<double>::quiet_NaN();
    }
    record("fitresult_state_max_abs_difference", maxDifference);
  }
'''


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ninja_command(target):
    output = subprocess.check_output([str(ENV / 'bin/ninja'), '-t', 'commands', target], cwd=BUILD, text=True)
    command = shlex.split(output.strip().splitlines()[-1])
    if command[:2] == [':', '&&']:
        command = command[2:]
    if command[-2:] == ['&&', ':']:
        command = command[:-2]
    assert '&&' not in command and ';' not in command
    return command


def main():
    started = time.time()
    DEST.mkdir(parents=True, exist_ok=True)
    manifest_file = DEST / 'build_manifest.json'
    if manifest_file.exists():
        raise RuntimeError('A completed isolated runtime already exists; refusing overwrite')
    source = DEST / 'source'
    for name in ('src', 'interface'):
        shutil.copytree(ORIGINAL / name, source / name, dirs_exist_ok=True)
    original_text = (ORIGINAL / 'src/GoodnessOfFit.cc').read_text()
    text = original_text.replace('#include <memory>', '#include <memory>\n#include <map>\n#include <limits>\n#include "../interface/ProfilingTools.h"')
    assert text.count('  minimn.minimize(verbose-2);') == 1
    text = text.replace('  minimn.minimize(verbose-2);', '  const bool nominalMinimizeOk = minimn.minimize(verbose-2);')
    assert text.count('  double nll_nominal = nominal_nll->getVal();') == 1
    text = text.replace('  double nll_nominal = nominal_nll->getVal();', '  double nll_nominal = nominal_nll->getVal();\n' + CAPTURE)
    (source / 'src/GoodnessOfFit.cc').write_text(text)
    (DEST / 'GoodnessOfFit.diff').write_text(''.join(difflib.unified_diff(original_text.splitlines(True), text.splitlines(True), fromfile='original/src/GoodnessOfFit.cc', tofile='diagnostic/src/GoodnessOfFit.cc')))
    copied_objects = {}
    for obj in (BUILD / 'CMakeFiles').rglob('*.o'):
        relative = obj.relative_to(BUILD)
        target = DEST / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(obj, target)
        copied_objects[str(relative)] = sha(obj)
    (DEST / 'lib').mkdir(exist_ok=True)
    (DEST / 'bin').mkdir(exist_ok=True)
    for name in ('libHiggsAnalysisCombinedLimit.rootmap', 'libHiggsAnalysisCombinedLimit_rdict.pcm'):
        shutil.copy2(BUILD / 'lib' / name, DEST / 'lib' / name)
    target = 'CMakeFiles/HiggsAnalysisCombinedLimit.dir/src/GoodnessOfFit.cc.o'
    compile_command = ninja_command(target)
    compile_command = [str(source / 'src/GoodnessOfFit.cc') if item == str(ORIGINAL / 'src/GoodnessOfFit.cc') else item for item in compile_command]
    link_command = ninja_command('lib/libHiggsAnalysisCombinedLimit.dylib')
    commands = [compile_command, link_command]
    live_files = [ORIGINAL / 'src/GoodnessOfFit.cc', ENV / 'bin/combine', ENV / 'lib/libHiggsAnalysisCombinedLimit.dylib', BUILD / 'lib/libHiggsAnalysisCombinedLimit.dylib']
    live_before = {str(p):sha(p) for p in live_files}
    with (DEST / 'build.log').open('w') as log:
        for command in commands:
            log.write(shlex.join(command) + '\n'); log.flush()
            subprocess.run(command, cwd=DEST, stdout=log, stderr=subprocess.STDOUT, check=True)
    # Copy the installed executable, then alter only its private dependency entry.
    shutil.copy2(ENV / 'bin/combine', DEST / 'bin/combine')
    private_library = DEST / 'lib/libHiggsAnalysisCombinedLimit.dylib'
    dependency_command = ['/usr/bin/install_name_tool', '-change', '@rpath/libHiggsAnalysisCombinedLimit.dylib', str(private_library), str(DEST / 'bin/combine')]
    subprocess.run(dependency_command, check=True)
    subprocess.run(['/usr/bin/codesign', '--force', '--sign', '-', str(DEST / 'bin/combine')], check=True)
    # Keep dependency environment setup exactly as the existing project runtime.
    wrapper = '#!/bin/zsh\nexec ' + shlex.join([str(ORIGINAL.parent / 'tools/pixi'), 'run', '--frozen', '--manifest-path', str(ORIGINAL / 'pixi.toml'), str(DEST / 'bin/combine')]) + ' "$@"\n'
    (DEST / 'combine').write_text(wrapper)
    (DEST / 'combine').chmod(0o755)
    live_after = {str(p):sha(p) for p in live_files}
    assert live_before == live_after
    result = dict(runtime_flag='GOF_TRACK_NOMINAL', native_sources_unchanged=True, live_before=live_before, live_after=live_after,
                  source_sha256=sha(source / 'src/GoodnessOfFit.cc'), diff_sha256=sha(DEST / 'GoodnessOfFit.diff'),
                  executable_sha256=sha(DEST / 'bin/combine'), library_sha256=sha(private_library),
                  reused_object_sha256=copied_objects, commands=commands + [dependency_command], seconds=time.time()-started,
                  wrapper=str(DEST / 'combine'), branches='gof_nominal_*; all Double_t, captured before saturated minimization',
                  no_additional_fits=True, no_model_changes=True)
    manifest_file.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('commands','reused_object_sha256','live_before','live_after')}, indent=2))


if __name__ == '__main__':
    main()
