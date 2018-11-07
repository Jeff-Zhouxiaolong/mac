function emd_addtag(obj,str)
//  add a tag to an object
// Calling Sequence
// emd_addtag(OBJ_HANDLE,STR)
//  Description
// Adds the tag STR to the object referrenced by OBJ_HANDLE
// When OBJ_HANDLE is an array of handles, STR is added to all corresponding
// objects.
//
// Rem: In order for this to work properly, the object's tag field must be a string
// containing keywords (or tags) separated by commas.
//
//
// See also
//  emd_hastag
//  emd_findtag
//  emd_rmtag
// Authors
// H. Nahrstaedt
// G.Rilling 12/2006 gabriel.rilling@ens-lyon.fr

inds = ~emd_hastag(obj,str);
if typeof(obj)=="handle"
  if(inds)
     if(isempty(obj.user_data))
       obj.user_data=str;
     else
       obj.user_data=[obj.user_data;str];
     end
  end
elseif typeof(obj)=="list"
  for l=1:length(obj)
     if(inds(l))
      if(isempty(obj(l).user_data))
	obj(l).user_data=str;
      else
	obj(l).user_data=[obj(l).user_data;str];
      end
    end
  end
end
// tags = get(obj(inds),'Tag');
//
// inds = find(inds);
// new_tags = cellfun(@(x)[x,',',str],tags,'UniformOutput',false);
// new_tags = cellfun(@(x)regexprep(x,'^,',''),new_tags,'UniformOutput',false);
// arrayfun(@(ind)set(obj(inds(ind)),'Tag',new_tags{ind}),1:length(inds));
endfunction
