#!/bin/tcsh

# give group permission to directories
find ./ -type d | xargs chmod g+s

# change group to 'gluex' for everything in path
chgrp -R gluex ./

# give group read permission
chmod g+rx -R ./
