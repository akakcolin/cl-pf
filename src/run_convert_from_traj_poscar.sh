for i in `cat filelist`
do
    cd $i
    for j in {0..19} 
    do
        cd $j
        ~/bin/convert_conf.x OUT.H5
        atomsk OUT.DB.xsf POSCAR
        cd ..
    done
    cd ..
done
