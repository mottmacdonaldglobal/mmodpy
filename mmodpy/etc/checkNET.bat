net session>nul 2>&1

:main
reg query "HKLM\SOFTWARE\Microsoft\NET Framework Setup\NDP\v4\Full"| findstr Version
exit