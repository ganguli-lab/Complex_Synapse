Get-CimInstance Win32_Process -Filter "name = 'python.exe'" | select ProcessId, CommandLine | Format-List
