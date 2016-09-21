% Now, the data from the test is ready for inspection
plot(tm, real(truth(:,k)), ';R(truth);',
	tm, imag(truth(:,k)), ';I(truth);',
	tm, real(output(:,k)), ';R(Out);', tm, imag(output(:,k)), ';I(Out);');
	grid on;

