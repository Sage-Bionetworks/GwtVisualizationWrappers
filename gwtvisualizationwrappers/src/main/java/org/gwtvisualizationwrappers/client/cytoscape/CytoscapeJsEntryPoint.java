package org.gwtvisualizationwrappers.client.cytoscape;

/*
 * #%L
 * GwtVisualizationWrapper
 * %%
 * Copyright (C) 2015 - 2016 GwtVisualizationWrapper
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

import com.google.gwt.core.client.EntryPoint;
import com.google.gwt.core.client.ScriptInjector;

/**
 * Provides script injection for cytoscape
 * 
 */
public class CytoscapeJsEntryPoint implements EntryPoint {


    /**
     * Check to see if jQuery is loaded already
     *
     * @return true is jQuery is loaded, false otherwise
     */
    private native boolean isjQueryLoaded() /*-{
        return (typeof $wnd['jQuery'] !== 'undefined');
    }-*/;
    
    /**
     * Check to see if Utility JS function is loaded already.
     * 
     */
    private native boolean isCytoscapeUtilsLoaded() /*-{
        return typeof $wnd['jQuery'].fn.createPlainObject !== 'undefined'
    }-*/;

    /** {@inheritDoc} */
    @Override
    public void onModuleLoad() {
        if (!isjQueryLoaded()) {
            ScriptInjector.fromString(CytoscapeJsClientBundle.INSTANCE.jQuery().getText())
                    .setWindow(ScriptInjector.TOP_WINDOW)
                    .inject();
        }
        if (!isCytoscapeUtilsLoaded()) {
        	ScriptInjector.fromString(CytoscapeJsClientBundle.INSTANCE.sageCytoscapeUtils().getText())
	            .setWindow(ScriptInjector.TOP_WINDOW)
	            .inject();
        }

	    

    }
    
}
