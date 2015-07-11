GwtVisualizationWrappers contain GWT wrappers for common js visualization libraries.  
Plan to include one gwt module for [Cytoscape JS](http://js.cytoscape.org/)
and
a genome browser.


##How to use
In GWT, after a container widget (div) has been attached to the DOM, we can start to render the visualization in that element.

		containerWidget.addAttachHandler(new AttachEvent.Handler() {
			@Override
			public void onAttachOrDetach(AttachEvent event) {
				if (event.isAttached()) {
					//render visualization
					String id = "cy1";
					containerWidget.getElement().setId(id);
					String cytoscapeGraphJson = "{\"elements\":{\"nodes\":[{\"data\":{\"id\":\"foo\"},\"position\":{\"x\":120,\"y\":120}},{\"data\":{\"id\":\"bar\"},\"position\":{\"x\":110,\"y\":110}},{\"data\":{\"weight\":100},\"group\":\"nodes\",\"position\":{\"x\":100,\"y\":100},\"selected\":true,\"selectable\":true,\"locked\":true,\"grabbable\":true}],\"edges\":[{\"data\":{\"id\":\"baz\",\"source\":\"foo\",\"target\":\"bar\"}}]},\"style\":[{\"selector\":\"node\",\"style\":{\"content\":\"data(id)\"}}]}";
					new CytoscapeGraph242().show(id,  cytoscapeGraphJson);
				};
			}
		});


##How to add a new js library to this project
We need to version each api call, to allow Sage Bionetworks to update the visualization version without breaking all references.
To support this, please run the external js through a compiler to obfuscate, and version each visualization call.
Also, script injection should be done lazily (when visualization is invoked).

1.  Download closure compiler from Google:
https://developers.google.com/closure/compiler/docs/gettingstarted_app
2.  Download external js library.
3.  Create new js library file.  For example:
java -jar compiler.jar --language_in ECMASCRIPT5_STRICT cytoscape.js --js_output_file 2.4.2-cytoscape.closurecompiled.js
4.  Add js as a TextResource (see CytoscapeJsClientBundle)
5.  Create the EntryPoint to the module (that loads any small shared resources up front).
6.  Create versioned visualization java object, that lazily loads the correct js resource and deals with version collisions.  See CytoscapeGraph242. 