LIBNAME=logspline
LD=ld

OBJS= logsplin.o mpack.o

$(LIBNAME).so:	$(OBJS)
	@$(LD) $(SHLIBLDFLAGS) -o $(LIBNAME).so $(OBJS)

clean:
	@rm -f *.o *.so

realclean: clean
