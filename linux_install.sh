#!/bin/bash
set -e
if ! python_loc="$(type -p python3)" || [[ -z $python_loc ]]; then
	echo "'python3' command not found, Python version >= 3.9 and < 3.12 is required."
fi
if ! git_loc="$(type -p git)" || [[ -z $git_loc ]]; then
	echo "'git' command not found, please install git and make it available on PATH to install."
fi
read -p "Do you wish to install to a specific folder? (Press Enter to use the default current directory) " install_folder
if [[ -z $install_folder ]]; then
	install_folder=$(pwd)
else
	if [ "$install_folder" == "." ]; then
		install_folder=$(pwd)
	elif [ "$install_folder" == ".."]; then
		install_folder=$(dirname $(pwd))
	else
		install_folder=$(cd $(dirname $install_folder); pwd)/$(basename $install_folder)
fi
mkdir -p $install_folder/flower_pot
echo "Cloning git repository..."
git clone --depth 1 --single-branch https://github.com/toboooo/flower_pot.git $install_folder/flower_pot
cd $install_folder
echo "Creating python environment..."
python3 -m venv flwpt
source flwpt/bin/activate
echo "Installing required packages..."
pip install --no-cache-dir -r flower_pot/requirements.txt
echo "Creating executable file in ${HOME}/Desktop folder..."
cat > $HOME/Desktop/flowerpot.sh <<!
#!/bin/bash
cd ${install_folder}
source flwpt/bin/activate
cd flower_pot
python3 flowerpot.py
!
chmod +x $HOME/Desktop/flower_pot.sh
