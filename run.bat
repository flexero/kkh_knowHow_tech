@echo off
rem Runs the engine offline. Order: portable Python in ..\_runtime\python, then py launcher, then python.
setlocal
set "HERE=%~dp0"
set "ENGINE=%HERE%hvac_engine.py"
set "RT=%HERE%..\_runtime\python\python.exe"
if exist "%RT%" ( "%RT%" "%ENGINE%" %* & goto end )
where py >nul 2>nul && ( py -3 "%ENGINE%" %* & goto end )
where python >nul 2>nul && ( python "%ENGINE%" %* & goto end )
echo Python 3.6+ not found. Install Python, or unzip a portable Python into _runtime\python\
:end
pause
