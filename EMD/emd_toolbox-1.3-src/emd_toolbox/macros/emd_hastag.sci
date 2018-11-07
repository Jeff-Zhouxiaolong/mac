function bool=emd_hastag(obj,str)
//  tests if an object has a specific tag
// Calling Sequence
// BOOL = emd_hastag(OBJ_HANDLE,STR)
// Description
// Tests if the object corresponding to OBJ_HANDLE has the tag STR.
// When OBJ_HANDLE is an array of handles emd_hastag returns a logical array of
// the same size.
//
// Rem: In order for this to work properly, the object's tag field must be a string
// containing keywords (or tags) separated by commas.
//
//
// See also
//  emd_addtag
//  emd_findtag
//  emd_rmtag
// Authors
// H. Nahrstaedt - Aug 2010
// G.Rilling 12/2006 gabriel.rilling@ens-lyon.fr

if (typeof(obj)=="handle")
tag = obj.user_data;
if isempty(tag)
  bool=%f;
elseif type(tag)==10
  bool = ~isempty(regexp(tag,['/'+str+'/'],'o'));
elseif type(tag)==15
  bool=%f;
  for k=1:length
   bool=or([bool, ~isempty(regexp(tag(k),['/'+str+'/'],'o'))]);
  end
end
elseif (typeof(obj)=="list")
  bool=ones(1,length(obj))*%f;
  for l=1:length(obj)
    tag = obj(l).user_data;
    if isempty(tag)
      bool(l)=%f;
    elseif type(tag)==10
      bool(l) = ~isempty(regexp(tag,['/'+str+'/'],'o'));
    elseif type(tag)==15
      bool(l)=%f;
      for k=1:length
      bool(l)=or([bool(l), ~isempty(regexp(tag(k),['/'+str+'/'],'o'))]);
      end
    end
  end
end
endfunction
