#!/bin/bash
BOOST_BRANCH="boost-1.70.0" git submodule foreach\
  'case $name in *"boost"*) git checkout $BOOST_BRANCH ;; *) ;; esac'