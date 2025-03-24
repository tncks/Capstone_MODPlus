package processedDB;

import java.util.Collections;
import java.util.LinkedList;

public class CandidateContainer {
	
	int size;
	final MODPeptide[] modlist;
	
	public int size() { return size; }
	public MODPeptide[] getList() { return modlist; }
	
	public CandidateContainer( LinkedList<TagPeptide> hmap, TagTrie trie ){ //multi-mod
		hmap.sort(new TagPeptComparator());

		size  = 0;
		modlist = new ChainTagPeptide[hmap.size()];

		while(!hmap.isEmpty()){
			TagPeptide parent = hmap.getFirst();
			ChainTagPeptide ctp = new ChainTagPeptide(parent.pStart, parent.pEnd, parent.mTag);
			ctp.setConservedRegion(parent.pLeft, parent.pRight);
			hmap.removeFirst();
			while(!hmap.isEmpty()){
				TagPeptide entry = hmap.getFirst();	
				if( !ctp.extend(parent, entry, trie) ){
					ctp.arrangeTags();
					modlist[size++] = ctp;
					ctp = new ChainTagPeptide(entry.pStart, entry.pEnd, entry.mTag);
					ctp.setConservedRegion(entry.pLeft, entry.pRight);
				}	
				parent = entry;
				hmap.removeFirst();
			}	
			ctp.arrangeTags();
			modlist[size++] = ctp;
		}		
	}
	
	public CandidateContainer( LinkedList<MODPeptide> hmap ){ //one-mod
		Collections.sort( hmap );
	
		size  = 0;
		modlist = new MODPeptide[hmap.size()];
		
		while(!hmap.isEmpty()){
			
			MODPeptide parent = hmap.getFirst();
			hmap.removeFirst();
			while(!hmap.isEmpty()){
				MODPeptide entry = hmap.getFirst();			
				if( !parent.extend(entry) ){
					modlist[size++] = parent;				
					parent = entry;
				}			
				hmap.removeFirst();
			}	
			modlist[size++] = parent;
		}		
	}
}
