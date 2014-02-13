@echo off

for /F "delims=" %%f in ('dir /s /b *.svg') do (
"C:\Program Files (x86)\Inkscape\inkscape" -f %%~nf.svg -d 300 -e %%~nf.png
)