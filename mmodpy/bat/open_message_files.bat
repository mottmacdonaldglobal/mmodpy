net session>nul 2>&1

:main
for /F "delims=" %%a in ('dir /s /b /a-d "mes0000"') DO "C:\Program Files (x86)\Notepad++\notepad++.exe" "%%a"
exit