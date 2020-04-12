bin_SCRIPTS = gist-rfe gist-sigmoid   \
gist2html  gist-score-svm

CLEANFILES = ${bin_SCRIPTS}

EXTRA_DIST = gist-rfe.in gist-sigmoid.in gist2html.in gist-score-svm.in

gist-rfe: 
	${PERL} -c gist-rfe.in
	rm -f gist-rfe
	echo "#!" ${PERL} > gist-rfe
	cat gist-rfe.in >> gist-rfe
	chmod a+x gist-rfe

gist-score-svm: 
	${PERL} -c gist-score-svm.in
	rm -f gist-score-svm
	echo "#!" ${PERL} > gist-score-svm
	cat gist-score-svm.in >> gist-score-svm
	chmod a+x gist-score-svm

gist2html: 
	${PERL} -c gist2html.in
	rm -f gist2html
	echo "#!" ${PERL} > gist2html
	cat gist2html.in >> gist2html
	chmod a+x gist2html

gist-sigmoid: 
	${PERL} -c gist-sigmoid.in
	rm -f gist-sigmoid
	echo "#!" ${PERL} > gist-sigmoid
	cat gist-sigmoid.in >> gist-sigmoid
	chmod a+x gist-sigmoid

