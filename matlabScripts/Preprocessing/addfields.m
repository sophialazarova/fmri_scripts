function updated_s = addfields(s, s_fields)

if ~isstruct(s_fields); error('s_fields must be a struct'); end
updated_s = s;

for namefield = fieldnames(s_fields)'
    updated_s.(namefield{:}) = s_fields.(namefield{:});
end