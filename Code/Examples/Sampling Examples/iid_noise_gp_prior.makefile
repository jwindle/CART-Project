# makefile.iid_noise_gp_prior
# To run: make -f makefile.iid_noise_gp_prior

iid_noise_gp_prior.out : iid_noise_gp_prior.o
	g++ -static iid_noise_gp_prior.o -lgsl -lgslcblas -lm \
	-o iid_noise_gp_prior.out

iid_noise_gp_prior.o : iid_noise_gp_prior.c
	g++ -Wall -c iid_noise_gp_prior.c

clean:
	rm iid_noise_gp_prior.o iid_noise_gp_prior.out