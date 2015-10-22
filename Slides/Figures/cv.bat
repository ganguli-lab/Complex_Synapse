@echo off

for /F "delims=" %%f in ('dir /s /b *-generated.pdf') do (
pdfpagegroup %%~nf.pdf
)