ACLOCAL_AMFLAGS = -I config
SUBDIRS=include src test doxygen examples
EXTRA_DIST = doxy.conf.in

########## Docs

.PHONY: doc doc-dev

if DOXY
doc: doxy.conf
	doxygen doxy.conf

doc-dev: doxy.conf
	(cat doxy.conf; \
	echo "ENABLED_SECTIONS=DEV"; \
	echo "OUTPUT_DIRECTORY=doc-dev"; \
	echo "PREDEFINED = FAAST_TIMINGS"; \
	echo "SHOW_USED_FILES = YES"; \
	echo "SHOW_FILES = YES") | doxygen -

dist-doc: doc
	tar czf $(PACKAGE)-$(VERSION)-refman.tgz --transform='s|^doc|$(PACKAGE)-$(VERSION)-refman|' doc/html

else
doc doc-dev dist-doc:
	@echo "Install doxygen to generate the documentation."
endif
