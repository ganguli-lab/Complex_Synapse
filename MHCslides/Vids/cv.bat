@echo off

for /F "delims=" %%f in ('dir /s /b *.svg') do (
"C:\Program Files (x86)\Inkscape\inkscape" -f %%~nf.svg -w 1600 -h 736 -e %%~nf.png
)