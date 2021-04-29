net session>nul 2>&1

:main
for /f %%F in ('dir /b /a-d ^| findstr /vile ".key .sh machines .bat"') do del "%%F"
exit