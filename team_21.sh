#!/bin/bash
path_to_train=$1
path_to_test=$2
path_to_store=$3

matlab  -nosplash  -nodesktop  -nodisplay  -nojvm  -r  "Directory('$path_to_train', '$path_to_test','$path_to_store');quit;"

