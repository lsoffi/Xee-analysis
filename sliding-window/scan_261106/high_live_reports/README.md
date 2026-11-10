# Live partial high-SR validation plots

Run `python3 -u watch.py` on lxplus with numpy, scipy and matplotlib.
The three independent workers read the existing nominal campaigns; no fits or
Condor jobs are submitted. Each successful report is refreshed after 300 seconds.
Reports are published to `www/Xee/261106/scan_SR_CR/high_live/` with CSV data,
counts and UTC timestamps. Incomplete samples are provisional, including possible
completion-time selection effects. The bias report retains the existing branch
validity and boundary criteria. Coverage reads saved finite CLs outputs, not a
new independent minimizer audit. GoF uses the SR-only statistic from simultaneous
SR+CR fits. No nuisance or model settings are changed.

For resilience to AFS, copy these scripts to a dedicated directory on the login
host's local disk before starting with `setsid nohup`. The local lock prevents a
second updater on the same host. `watch_status.json` records host, PID, worker
status, and publication time. The web heartbeat must advance; an alive PID alone
is not proof of progress. The watcher can survive an SSH disconnect, but not a
host reboot or login-node cleanup. It does not send chat notifications.

Publication retains the previous plots when a report fails. Inspect each worker
log before repairing; do not loosen validity criteria to produce a figure.
