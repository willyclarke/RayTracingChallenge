set makeprg=ninja\ -C\ $PWD/../build/RayTracingChallenge
let g:chromatica#dotclangfile_search_path = $PWD

" path to directory where libclang.so can be found
let g:ncm2_pyclang#library_path = '/usr/local/opt/llvm/lib'
" or path to the libclang.so file
"let g:ncm2_pyclang#library_path = '/usr/lib64/libclang.so.5.0'

" a list of relative paths for compile_commands.json
let g:ncm2_pyclang#database_path = [
            \ 'compile_commands.json',
            \ 'build/compile_commands.json',
            \ '../build/RayTracingChallenge/compile_commands.json'
            \ ]
" a list of relative paths looking for .clang_complete
let g:ncm2_pyclang#args_file_path = ['.clang_complete']


