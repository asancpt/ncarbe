roxygen:
	Rscript -e "roxygen2::roxygenise()"

pkgdown:
	Rscript -e "pkgdown::build_site()"

readme:
	Rscript -e "rmarkdown::render('README.Rmd', output_format = 'github_document')"
