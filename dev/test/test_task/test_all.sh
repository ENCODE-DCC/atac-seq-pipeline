#!/bin/bash
set -e

for wdl in test_*.wdl
do
  json=${wdl%.*}.json
  result=${wdl%.*}.result.json
  echo "./test.sh ${wdl} ${json} ${1}"
  #./test.sh ${wdl} ${json} ${1}
  #python -c "import sys; import json; data=json.loads(sys.stdin.read()); sys.exit(int(not data[u'match_overall']))" < ${result}
  #rm -f ${result}
done
