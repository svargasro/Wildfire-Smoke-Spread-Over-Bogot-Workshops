GCC = g++
CFLAGS = -O3 --std=c++17


TARGET = LB_AdvectionDiffusion_D2Q9
# TARGET = $(patsubst %.cpp,%,$(wildcard *.py))

all: ejecutar gif 

$(TARGET).x: $(TARGET).cpp
	$(GCC) $(CFLAGS) -o $@ $<

ejecutar: $(TARGET).x
	./$< 0.016	

gif:
#	convert -delay 25 -loop 0 $(shell ls frames/*.png | sort -V) AdvectionDifusion.gif
	$(eval FRAMES=$(shell ls frames/*.png | sort -V))
	$(eval COUNT=$(shell echo "$$(ls frames/*.png | wc -l)"))
	$(eval BLOCK_SIZE=45)
	$(eval GIF_COUNT=$(shell echo $$(($(COUNT)/$(BLOCK_SIZE)))))

	for i in $$(seq 0 $$(($(GIF_COUNT)-1))); do \
		START=$$(($$(($$i*$(BLOCK_SIZE)))+1)); \
		END=$$(($$START+$(BLOCK_SIZE)-1)); \
		convert -delay 10 -loop 0 $$(ls frames/*.png | sort -V | sed -n "$$START,$$END p") AdvectionDifusion_$$i.gif; \
	done

.PHONY: clean

clean:
	rm *.x  *.gp
	rm -r frames data