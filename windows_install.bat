@echo off
python --version >nul 2>nul
if not %ERRORLEVEL% equ 0 echo "'python' command not found, Python version >= 3.9 and < 3.12 is required." && exit /b
git --version >nul 2>nul
if not %ERRORLEVEL% equ 0 echo "'git' command not found, please install git and make it available on %PATH% to install." && exit /b
set /p "install_folder=Do you wish to install to a specific folder? (Press Enter to use the default current directory) "
if %install_folder% == %^install_folder% (
	set install_folder=%cd%
)
else (
	set "install_folder=%install_folder:/=\%"
	pushd %install_folder%
	set install_folder=%cd%
	popd
)
if not %ERRORLEVEL% equ 0 exit /b
mkdir %install_folder%\flower_pot
if not %ERRORLEVEL% equ 0 exit /b
echo "Cloning git repository..."
git clone --depth 1 --single-branch https://github.com/toboooo/flower_pot.git %install_folder%\flower_pot
if not %ERRORLEVEL% equ 0 exit /b
cd %install_folder%
echo "Creating python environment..."
python -m venv flwpt
if not %ERRORLEVEL% equ 0 exit /b
call flwpt\Scripts\activate
if not %ERRORLEVEL% equ 0 exit /b
echo "Installing required packages..."
pip install --no-cache-dir -r flower_pot\requirements.txt
if not %ERRORLEVEL% equ 0 echo "Problem installing Flower Pot's dependencies. Is the Python version 3.9 or higher and 3.11 or lower?" && exit /b
echo "Creating executable file in Desktop folder (%HOMEDRIVE%%HOMEPATH%\Desktop)..."
set "desktop_file=%HOMEDRIVE%%HOMEPATH%\Desktop\flowerpot.bat"
echo @echo off>%desktop_file%
if not %ERRORLEVEL% equ 0 exit /b
echo cd %install_folder%>>%desktop_file%
echo call flwpt\Scripts\activate>>%desktop_file%
echo cd flower_pot>>%desktop_file%
echo python flowerpot.py>>%desktop_file%
