#!/bin/bash

# Copy Python scripts to the bin directory
cp $SRC_DIR/scripts/puppy-align $PREFIX/bin/puppy-align
cp $SRC_DIR/scripts/puppy-primers $PREFIX/bin/puppy-primers
cp $SRC_DIR/scripts/puppy-GUI $PREFIX/bin/puppy-GUI

# Make them executable
chmod +x $PREFIX/bin/puppy-align
chmod +x $PREFIX/bin/puppy-primers
chmod +x $PREFIX/bin/puppy-GUI
