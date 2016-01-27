
export PATH=$PATH:/usr/bin/:/bin/

if [ "$?" != "0" ]; then
	echo "Error: Build failed"
	exit 1
fi

PACKAGE=mdy

find $PACKAGE -type d -name __pycache__ -exec rm -rf {} \; -print || true

echo "Installing into $PREFIX"

for PYTHON_VER in 3.4 3.5; do
	DIR="$SP_DIR/../../python${PYTHON_VER}/site-packages"

	if [ "$DIR" != "" ]; then
		mkdir -p "$DIR"
	fi
	if [ -e "$DIR" ]; then
		pwd
		ls
		cp -r python/$PACKAGE  $DIR/
    rm -rf $DIR/$PACKAGE/data
	  rm -rf $DIR/$PACKAGE/.idea



	  for D in $PACKAGE ; do
	  	cd $DIR/$D
	 		TMP=$(mktemp /tmp/XXXXXX) 
			for T in $(find . -name "*.py"); do
            echo "#!/usr/bin/env python" > $TMP
  	  			echo "# (c) 2015 Acellera Ltd" >> $TMP
    				echo "# All Rights Reserved" >> $TMP
	    			echo "# Distributed under HTMD Software Academic License Agreement v1.1" >> $TMP
			  	  echo "# No redistribution in whole or part" >> $TMP 
						echo "#" >> $TMP
						cat $T >> $TMP
						cp $TMP $T
			done
			rm $TMP
      cd -
		done
	else
		echo "Error: SP_DIR not defined"
		exit 1
	fi
	chmod -R a+rX $DIR
done


exit 0
