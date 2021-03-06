########################################################################
# This Makefile can be used to execute your program once it is compiled.
# It can contain commands that generate images in your submitted 
# output directory from the ones in your input directory.
########################################################################



########################################################################
# Executable name
########################################################################

EXE=src/photonmap



########################################################################
# This command tells the make file how to execute src/lap
# to produce output from basic rules.
########################################################################

output/%.jpg: input/%.scn
	$(EXE) $< $@



########################################################################
# "make all" or simply "make" should make all output images
########################################################################

all: \
	$(EXE) \
        output/cornell.jpg \
        output/cornell2.jpg \
        output/cos526.jpg \
        output/dirlight1.jpg \
        output/dirlight2.jpg \
        output/reflection.jpg \
        output/specular.jpg \
        output/spotlight2.jpg 

clean:
	cd src; make clean
	rm -f output/*

$(EXE): src
	cd src; make



