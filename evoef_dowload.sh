path=`python -c 'import site; print(site.getsitepackages()[0])'`
path=${path}/EESMHM
echo $path
cd $path
wget https://github.com/tommyhuangthu/EvoEF/archive/refs/heads/master.zip
unzip master.zip
rm master.zip
cd EvoEF-master
chmod +x build.sh
./build.sh

