
project(
'rapidsnark', 
'cpp',
default_options: [
  'buildtype=release',
  'debug=true',
  'cpp_std=c++17'
  ]
)

# for ffiasm
#add_languages('nasm')

# b/c ffiasm uses these flags
# add_project_arguments('-DUSE_ASM', language : 'cpp')
add_project_arguments('-DARCH_ARM', language : 'cpp')


subdir('src')



