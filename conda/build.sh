#!/bin/bash

# Copy Python scripts to the bin directory
cp $SRC_DIR/scripts/puppy-align.py $PREFIX/bin/puppy-align
cp $SRC_DIR/scripts/puppy-primers.py $PREFIX/bin/puppy-primers

# Copy additional scripts
cp $SRC_DIR/scripts/MMseqs2.sh $PREFIX/bin/MMseqs2
cp $SRC_DIR/scripts/UniqueGenesPlotting.R $PREFIX/bin/UniqueGenesPlotting

# Make them executable
chmod +x $PREFIX/bin/puppy-align
chmod +x $PREFIX/bin/puppy-primers
chmod +x $PREFIX/bin/MMseqs2
chmod +x $PREFIX/bin/UniqueGenesPlotting