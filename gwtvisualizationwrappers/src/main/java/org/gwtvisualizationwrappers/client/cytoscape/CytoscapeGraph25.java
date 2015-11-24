package org.gwtvisualizationwrappers.client.cytoscape;

import com.google.gwt.core.client.ScriptInjector;

public class CytoscapeGraph25 {
    
	/**
     * Check to see if Cytoscape JS version has been loaded already.
     * 
     * @return true if Cytoscape is loaded, false otherwise.
     */
    private native boolean isCytoscape25Loaded() /*-{
        return typeof $wnd['jQuery'].fn.cytoscape25 !== 'undefined'
    }-*/;
    
	private static native void _init25() /*-{
		$wnd.cytoscape25 = $wnd.cytoscape;
		$wnd.cytoscape = undefined;
	}-*/;
	
	private static native void _initGraph(String containerId, String cyjs, String styleJson) /*-{
		var containerElement = $doc.getElementById(containerId);
		function readyFunction() {
			console.log('Cytoscape graph ready');
		}
	
		var options = $wnd.createPlainObject(cyjs, styleJson, containerElement, readyFunction);
		
		$wnd.cytoscape = $wnd.cytoscape25;
		$wnd.cytoscape(options); 
		$wnd.cytoscape = undefined;
	}-*/;
	
	/**
	 * Construct and show a cytoscape graph.
	 * 
	 * @param containerId Element ID to put the graph into.
	 * @param cytoscapeGraphJson Exported JSON from Cytoscape.
	 * http://wiki.cytoscape.org/Cytoscape_3/UserManual#Cytoscape_3.2BAC8-UserManual.2BAC8-CytoscapeJs.Data_Exchange_between_Cytoscape_and_Cytoscape.js

{
	elements:{
		nodes:[],
		edges:[]
	}
	style:[{
			selector: 'node',
			style: {
			}
	}]
}
	 */
	public void show(String containerId, String cyjs, String styleJson) {
		//lazy load the cytoscape.js source
		if (!isCytoscape25Loaded()) {
		    ScriptInjector.fromString(CytoscapeJsClientBundle.INSTANCE.cytoscape2_5().getText())
		        .setWindow(ScriptInjector.TOP_WINDOW)
		        .inject();
		    _init25();
		}

		_initGraph(containerId, cyjs, styleJson);
	}



}
