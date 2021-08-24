@echo off
@REM python.exe is attempted first; this is because, on Windows, python.exe is the "managed" install,
@REM   and as such is the most likely candidate for the correct Python installation.
@REM If the managed executable (python.exe) is not v3, then we can assume that python3.exe is the right one.
for /f "tokens=1" %%i in ('powershell.exe ^(python.exe --version^).Split^(^)[1][0]') do set _PY_VERSION_=%%i
if %_PY_VERSION_% == 3 (
    python.exe %~dp0\hybrid_c_wrapper.py %*
) else (
    python3.exe %~dp0\hybrid_c_wrapper.py %*
)
