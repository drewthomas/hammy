hammy: hammy.c nwgks.[ch] data.[ch] Makefile
	gcc -Os -s -Wall -Wextra -ansi -pedantic -std=c99 *.c -lm -o hammy

nwgks: nwgks.[ch] data.[ch] Makefile
	gcc -Os -s -Wall -Wextra -ansi -pedantic -std=c99 \
		-DMAKE_STANDALONE_NWGKS nwgks.c data.c -lm -o nwgks

.PHONY nwgksloo_output \
output/nwgksloo_JS.dat output/nwgksloo_CO2.dat output/nwgksloo_RVN63.dat \
output/nwgksloo_RPDI.dat output/nwgksloo_diffract.dat \
output/nwgksloo_cosine.dat output/nwgksloo_multicos.dat \
output/nwgksloo_zigzag.dat: nwgksloo data/*.dat
	./nwgksloo < data/JS_recount_money.dat > output/nwgksloo_JS.dat
	./nwgksloo < data/NOAA_ESRL_CO2.dat > output/nwgksloo_CO2.dat
	./nwgksloo < data/RVN63_mil_bivar.dat > output/nwgksloo_RVN63.dat
	./nwgksloo < data/UK_RPDI.dat > output/nwgksloo_RPDI.dat
	./nwgksloo < data/diffraction.dat > output/nwgksloo_diffract.dat
	./nwgksloo < data/synthetic/cosine.dat > output/nwgksloo_cosine.dat
	./nwgksloo < data/synthetic/zigzag.dat > output/nwgksloo_zigzag.dat
	./nwgksloo < data/synthetic/multi_cosine.dat > output/nwgksloo_multicos.dat
