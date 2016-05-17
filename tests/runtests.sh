#!/bin/bash

echo -e "\033[95mRunning unit tests:\033[0m"

for i in bin/test_*
do
	if test -f $i
	then
		if $VALGRIND ./$i 2>> ./tests.log
		then
			echo $i PASS
		else
			echo "ERROR in test $i: here's ./tests.log"
			echo "------"
			tail ./tests.log
			exit 1
		fi
		if [ ! -z "$VALGRIND" ]
		then
			for log in `ls /tmp/valgrind-*.log`;
			do
				cmdstr=$(head -n 4 $log | tail -n 1)
				if [[ $cmdstr == *"$i"* ]]
				then
					tailstr=$(tail -n 1 $log)
					if [[ $tailstr == *"0 errors from 0 contexts"* ]]
					then
						echo -e "\033[92mVALGRIND GOOD\033[0m"
					else
						echo -e "\033[91mVALGRIND BAD\033[0m"
					fi
				fi
			done
		fi
	fi
done

echo ""
