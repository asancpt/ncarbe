roxygen:
	Rscript -e "roxygen2::roxygenise()"

pkgdown:
	Rscript -e "pkgdown::build_site()"

github:
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
