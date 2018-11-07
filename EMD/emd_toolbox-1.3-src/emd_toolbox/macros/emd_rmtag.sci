function emd_rmtag(obj,str)
//  removes tag from object
// Calling Sequence
// emd_rmtag(OBJ_HANDLE,STR)
// Description
// Removes the tag STR from the object referrenced by OBJ_HANDLE
// When OBJ_HANDLE is an array of handles, the tag is removed from all
// corresponding objects.
//
// Rem: In order for this to work properly, the object's tag field must be a string
// containing keywords (or tags) separated by commas.
//
// See also
//  emd_addtag
//  emd_hastag
//  emd_findtag
// Authors
// H. Nahrstaedt - Aug 2010
// G.Rilling 12/2006 gabriel.rilling@ens-lyon.fr


if or(~emd_hastag(obj,str))
    warning('emd_rmtag:warning','no such tag in object')
end
inds = ~emd_hastag(obj,str);
if typeof(obj)=="handle"
  if(inds)
     if (size(obj.user_data,'*')==1)
           obj.user_data==[];
     else
       for l=1:size(obj.user_data,'*')
           if ~isempty(regexp(tag(k),['/'+str+'/'],'o'))
              obj.user_data(l)==[];
           end
       end
     end;
  end
elseif typeof(obj)=="list"
  for l=1:length(obj)
    if (inds(l))
     if (size(obj(l).user_data,'*')==1)
           obj(l).user_data==[];
     else
       for l=1:size(obj(l).user_data,'*')
           if ~isempty(regexp(tag(k),['/'+str+'/'],'o'))
              obj(l).user_data(l)==[];
           end
       end
     end;
    end
  end
end



// arrayfun(@emd_rmtag1,1:length(obj));
//
//   function emd_rmtag1(ind)
//     tag = get(obj(ind),'Tag');
//     tag = regexprep(tag,['(\W|^)',str,'(\W|$)'],'$1$2');
//     tag = regexprep(tag,',,',',');
//     tag = regexprep(tag,'^,|,$','');
//     set(obj(ind),'Tag',tag);
//   end

endfunction
