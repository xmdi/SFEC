CC := g++
FLAGS := -Ofast -std=c++14
SRCDIR := src
BUILDDIR := bin
DATDIR := dat

core: dirs
	$(CC) $(SRCDIR)/main.cpp -o $(BUILDDIR)/core $(FLAGS)
test:  
	$(BUILDDIR)/core
	@echo Plotting Test 1, Rotate with Mouse, Enter to Continue
	@gnuplot -c $(SRCDIR)/view.gp 'testModel_1' -p
	@echo Plotting Test 2, Rotate with Mouse, Enter to Continue
	@gnuplot -c $(SRCDIR)/view.gp 'testModel_2' -p
	@echo Plotting Test 3, Rotate with Mouse, Enter to Continue
	@gnuplot -c $(SRCDIR)/view.gp 'testModel_3' -p

.PHONY: dirs clean
dirs:
	mkdir -p $(BUILDDIR) $(DATDIR)
clean: 
	rm -rf $(BUILDDIR) $(DATDIR) 
