# Restart state and why a clean run is different

`UOPT1/OPT_IEXEC` has six line-oriented values. Its **sixth (last) line** is
the restart flag read by `RI_EXEC` in `read.f`:

- `0` -> initialize `INCNUMBER=0`, `TTIME=0`, stress/plastic state and velocities
  from initial input; the first time integration increment is 1.
- `1` -> read `RESTART/PARA`, `STRESS_*`, `VELOCITY_*`, etc.; current example checkpoint begins at increment 35.

Note that `UOPT1/O_OEXEC`'s fields are different:

1. requested end time (`0.008` default experiment, original `0.6`)
2. Newton iteration limit (`220`)
3. initial dt (`2.5e-5`)
4–8. analysis/contact/relaxation/tolerance options (retained)

**Do not delete `UBIN1`, `UOPT1`, source files, or your historical checkpoints.**
Starting from zero is NOT achieved by merely deleting `RESTART` if the sixth
`OPT_IEXEC` record still equals 1.

Correct command: `bash scripts/run_from_zero_ucrt64.sh` on Windows UCRT64,
or `bash scripts/run_from_zero_linux.sh` on Linux/WSL.

Launchers move old `RESTART`, `OUT`, and log files into timestamped `run_backups/`,
create empty working directories, and explicitly set the flag 0 before executing.
You can pass `--end-time 0.6 --initial-dt 2.5e-5 --max-iterations 220`.

Check early output for `ADAPT_START   0    0.0` and first `ADAPT_STEP_START  1  0...`;
seeing `ADAPT_START 35` means you are *still restarting*.
