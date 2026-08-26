import nextflow.script.ChannelOut
import nextflow.Channel

class Checkpoint {

    /**
     * workflows emit a type ChannelOut, which is not a channel, 
     *   or a strict Map of Channels, but it has the fields that it can be kinda treated like that. 
     * This makes it into an actual Map<String, Channel>. 
     * except, that doesnt work either. 
     * Its only type 'Channel' at initialization time. At every time that matters it is
     *   some variation of DataflowSomething. So the value side of this map has a 
     *   cop-out type Object.  
     * Convert ChannelOut (or Map) into a real Map<String, Channel>.
     * - ChannelOut.channels is the named channels from emit:
     * - ChannelOut.target is the positional channels, if any. 
     *     Contains duplications of the named channels, but this removed those.
     *
     *  ArrayList is for cases where nextflow wants to expand a single ChannelOut 
     *   object into multiple workflow or function arguments. This causes argument count mismatches.
     *   so [workflow_name.out] prevents that, this can unwrap it.  
     */
    static Map<String, Object> normalizeOut(Object outObj) {
        if (outObj instanceof ChannelOut) {            
            Map<String,Object> m = new LinkedHashMap<>()
            if (outObj.channels) {
                m.putAll(outObj.channels as Map<String, Object>)
            }
            outObj.target?.eachWithIndex { ch, i ->
                String posKey = i.toString()
                // Is this positional channel identical to ANY named channel?
                boolean duplicate = m.any { namedKey, namedCh ->
                    // numeric positional key is synthetic only if the value is literally same object
                    !isNumeric(namedKey) && isSameChannel(namedCh, ch)
                }
                if (!duplicate) {
                    m[posKey] = ch
                }
            }
            return m
        }
        if (outObj instanceof Map) {
            return (Map<String, Object>) outObj
        }
        if (outObj instanceof ArrayList) {
            return normalizeOut(outObj[0])
        }
        throw new IllegalArgumentException("Expected ChannelOut or Map; got ${outObj?.getClass()}")
    }


    private static boolean isNumeric(Object key) {
        (key instanceof Number) || (key instanceof CharSequence && (key ==~ /^\d+$/))
    }

    private static boolean isSameChannel(Object a, Object b) {
        // Identity match (preferred) or equals match if identity unavailable
        a?.is(b) || a == b
    }

}

//

