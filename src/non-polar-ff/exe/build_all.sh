workdir=`pwd`

while read line           
do
    cd $workdir/$line
    sh compile.sh
    retcode=$?
    if [ "$retcode" -ne 0 ]; then
       echo "Error in $line."
       exit 1
    fi

done < dirlist.txt
