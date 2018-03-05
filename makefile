all:
	@mkdir -p dat
	@g++ main.cpp -o core -Ofast -std=c++14
test:
	@./core
	@echo Plotting Test 1, Rotate with Mouse, Enter to Continue
	@gnuplot -c ./view.gp 'testModel_1' -p
	@echo Plotting Test 2, Rotate with Mouse, Enter to Continue
	@gnuplot -c ./view.gp 'testModel_2' -p
	@echo Plotting Test 3, Rotate with Mouse, Enter to Continue
	@gnuplot -c ./view.gp 'testModel_3' -p
