package org.gwtvisualizationwrappers.client.cytoscape;
import com.google.gwt.dom.client.Element;
/*
 * #%L
 * GwtCytoscapeJs
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


public class CytoscapeGraph {
	/**
	 * Construct and show a cytoscape graph
	 * @param params
	 */
	public CytoscapeGraph(CytoscapeInitParams params) {
		show(params.container, params);
	}
	private native void show(Element e, CytoscapeInitParams params) /*-{
		//TODO: cytoscape(params) currently returns 'undefined' instead of the element.
		$wnd.jQuery(e).cytoscape(params);
	}-*/;
}
