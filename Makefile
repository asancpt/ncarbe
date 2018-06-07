publish: roxygen install clean readme pkgdown

roxygen:
	Rscript -e "roxygen2::roxygenise()"

pkgdown:
	Rscript -e "pkgdown::build_site()"

readme:
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"

install:
	cd ..; R CMD INSTALL --no-multiarch --with-keep.source ncarbe; cd ncarbe

clean:
	rm -rf docs
