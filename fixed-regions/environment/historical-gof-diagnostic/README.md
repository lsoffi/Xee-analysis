> **Archivio di provenienza.** Questo repository contiene solo diff e metadati: il binario e il wrapper `combine` citati sotto non sono distribuiti. Il setup corrente non applica questa patch. La descrizione seguente riguarda la macchina della campagna storica, non un runtime già disponibile nel checkout. Vedere anche `ARCHIVAL_NOTE.md`.

# Production-ready diagnostic runtime v2

This immutable version supersedes the first private build for production GoF diagnostics. Use the `combine` wrapper in this directory and add `--X-rtd GOF_TRACK_NOMINAL`.

The capture now includes **all ModelConfig nuisance parameters**, including `lumi_1` and `lumi_2`, as well as the previously captured local coefficients/yields, Υ scales and means, r, and family categories. All 118 branches use double precision. The source diff, commands, unchanged live-runtime hashes and private binary hashes are recorded in `build_manifest.json`.

The original toy007/m6.0 command was rerun into `validation` using the identical preserved workspace. Its native score and tracked r match the original exactly. Restoring the captured nominal state leaves no uncaptured floating PDF parameters or ModelConfig nuisances. The reconstructed exact extended deviance plus the native RooMultiPdf penalty agrees with the GoF score within 10⁻⁶. All input files and the live Combine runtime remain unchanged. `IMMUTABLE.json` records the final binary and validation hashes.

The first private runtime and its four score-parity checks remain preserved in the parent directory. It lacks two nuisance branches and must not be used for full-state reconstruction. The statistical calculation and diagnostic capture timing are otherwise the same; refer to the parent README for their description.

No additional minimization or Hesse call is introduced. GoF point validity uses convergence, EDM, positivity and bounds; covariance quality is recorded but may be 1 when no Hesse calculation was requested.
