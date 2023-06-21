#!/bin/bash

cd config_mainsel_UL16preVFP_ele
sframe_batch.py $1 parsedConfigFile_mainsel_UL16preVFP_ele.xml
cd ..
cd config_mainsel_UL16preVFP_muo
sframe_batch.py $1 parsedConfigFile_mainsel_UL16preVFP_muo.xml
cd ..

cd config_mainsel_UL16postVFP_ele
sframe_batch.py $1 parsedConfigFile_mainsel_UL16postVFP_ele.xml
cd ..
cd config_mainsel_UL16postVFP_muo
sframe_batch.py $1 parsedConfigFile_mainsel_UL16postVFP_muo.xml
cd ..

cd config_mainsel_UL17_ele
sframe_batch.py $1 parsedConfigFile_mainsel_UL17_ele.xml
cd ..
cd config_mainsel_UL17_muo
sframe_batch.py $1 parsedConfigFile_mainsel_UL17_muo.xml
cd ..

cd config_mainsel_UL18_ele
sframe_batch.py $1 parsedConfigFile_mainsel_UL18_ele.xml
cd ..
cd config_mainsel_UL18_muo
sframe_batch.py $1 parsedConfigFile_mainsel_UL18_muo.xml
cd ..
