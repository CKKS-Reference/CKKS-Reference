.code

mul_asm proc
	mov rax,rcx
	mul rdx
	mov [r8],rdx
	mov [r9],rax
	ret
mul_asm endp

div_asm proc
	mov rax,rdx
	mov rdx,rcx
	div r8
	mov [r9],rdx
	ret
div_asm endp

mod_asm proc
	mov rax,rdx
	mov rdx,rcx
	div r8
	mov rax,rdx
	ret
mod_asm endp

end