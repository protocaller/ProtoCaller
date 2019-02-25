#!/bin/bash

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    FILE="biosimspace_devel_latest_linux.run"
    if [[ ! -f $FILE ]]; then
        wget https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/ZH4wscDHe59T28yVJtrMH8uqifI_ih0NL5IyqxXQjSo/n/chryswoods/b/biosimspace_releases/o/biosimspace_devel_latest_linux.run
    fi
elif [[ "$OSTYPE" == "darwin"* ]]; then
    FILE = biosimspace_devel_latest_os.run
    if [[ ! -f $FILE ]]; then
        wget https://objectstorage.eu-frankfurt-1.oraclecloud.com/p/whcwfvWfndjA4RxupM-4gsVsjcdR0w5I9aP1RJKPruQ/n/chryswoods/b/biosimspace_releases/o/biosimspace_devel_latest_osx.run
    fi
else
    echo "Only Linux / MacOS systems are supported"
    exit 1
fi

chmod +x "$FILE"
{ echo ; } | "./$FILE" --nox11
source activate ~/biosimspace.app/
sleep 5
python setup.py -q install
