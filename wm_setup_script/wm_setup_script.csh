echo You are user: $USER

home_folder=home10

if ([ -d "/sciclone/$home_folder/$USER" ])
then
 	echo "Directory /sciclone/$home_folder/$USER exists."
else
	home_folder=home20
 	echo "Directory /sciclone/$home_folder/$USER exists."
fi

home_dir=/sciclone/$home_folder/$USER

echo $home_dir

#test_dir=wm_setup_script
test_dir=''

echo $home_dir/$test_dir

ln -s /sciclone/scr10 /$home_dir/$test_dir/scr
ln -s /sciclone/gluex10/$USER/ /$home_dir/$test_dir/gluex10
ln -s /local/scr/$USER /$home_dir/$test_dir/lscr
ln -s /sciclone/pscr/$USER /$home_dir/$test_dir/pscr

#ln -s /sciclone/scr10/ /sciclone/home20/wli08/


