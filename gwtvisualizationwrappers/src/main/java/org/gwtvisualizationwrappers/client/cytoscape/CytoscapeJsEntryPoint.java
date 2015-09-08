package org.gwtvisualizationwrappers.client.cytoscape;

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
