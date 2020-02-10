OBJDIR := objdir
OBJS := $(addprefix $(OBJDIR)/,foo.o bar.o baz.o)

$(OBJDIR)/%.o : %.f90
	@echo $@

all: $(OBJS) $(OBJDIR)/*.o 

$(OBJS): | $(OBJDIR)

$(OBJDIR):
	@mkdir $(OBJDIR)
