package org.gwtvisualizationwrappers.client.biodalliance;

import com.google.gwt.core.client.GWT;
import com.google.gwt.core.client.JavaScriptObject;
import com.google.gwt.core.client.JsonUtils;
import com.google.gwt.core.client.Scheduler;
import com.google.gwt.core.client.ScriptInjector;
import com.google.gwt.dom.client.Element;

/*
 * #%L
 * GwtVisualizationWrapper
 * %%
 * Copyright (C) 2015 GwtCytoscapeJs
 * %%
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *      http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * #L%
 */

public class Biodalliance013dev {
    
	/**
     * Check to see if Biodalliance JS version has been loaded already.
     * 
     * @return true if Biodalliance is loaded, false otherwise.
     */
    private native boolean isBiodalliance013DevLoaded() /*-{
    	return (typeof $wnd['Browser013Dev'] !== 'undefined');
    }-*/;
    
    
	private static native void _initBrowser(JavaScriptObject biodallianceBrowserConfig) /*-{
		$wnd.Browser = $wnd.Browser013Dev;
		var newBrowser = new $wnd.Browser(biodallianceBrowserConfig);
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
	public void show(JavaScriptObject biodallianceBrowserConfig) {
		//lazy load the cytoscape.js source
		if (!isBiodalliance013DevLoaded()) {
		    ScriptInjector.fromString(BiodallianceClientBundle.INSTANCE.biodalliance0_13dev().getText())
		        .setWindow(ScriptInjector.TOP_WINDOW)
		        .inject();
		    _init013Dev();
		}

		_initBrowser(biodallianceBrowserConfig);
	}

	private static native void _init013Dev() /*-{
		$wnd.Browser013Dev = $wnd.Browser;
	}-*/;


}
