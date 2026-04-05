-- .nvim.lua — Zig project local config
-- This file is trusted per-project; load with :set exrc

-- .nvim.lua (project root)
vim.api.nvim_create_autocmd("FileType", {
  pattern = "rust",
  callback = function()
    vim.bo.makeprg = "cargo test"
  end,
})

local ls = require("luasnip")
ls.cleanup() -- clears old snippets
local s = ls.snippet
local t = ls.text_node
local i = ls.insert_node

vim.keymap.set({ "i" }, "<C-Y>", function() ls.expand() end, { silent = true })
-- vim.keymap.set({ "i" }, "<C-K>", function() ls.expand() end, { silent = true })
-- vim.keymap.set({ "i", "s" }, "<C-L>", function() ls.jump(1) end, { silent = true })
vim.keymap.set({ "i", "s" }, "<C-J>", function() ls.jump(-1) end, { silent = true })

vim.keymap.set({ "i", "s" }, "<C-E>", function()
  if ls.choice_active() then
    ls.change_choice(1)
  end
end, { silent = true })

vim.keymap.set({ "i", "s" }, "<TAB>", function()
  if ls.expand_or_jumpable() then
    vim.fn.complete(0, {})  -- close completion menu
    ls.expand_or_jump()
  end
end, { silent = true })

ls.add_snippets("rust", {
  s("tst", {
    t({ "/// Chap x - ", "" }),
    t({ "#[test]", "fn " }),
    i(1, "test_chap_x_y"),
    t({ "() -> Result<(), String> {", "    " }),
    i(0),
    t({
      "let chk = Tuple::point(1.0, 1.0, 1.0).approx_eq(Tuple::point(1.0, 1.0, 1.0));",
      " "
    }),
    t({
      "if ",
    }),
    i(1, "chk"),
    t({
      " {",
      "    Ok(())",
      "} else {",
      "    Err(\"",
    }),
    i(2, ""),
    t({
      "\".into())",
      "}",
    }),
    t({ "", "}" }),
  }),
})

ls.add_snippets("rust", {
  s("rdoc", {
    t({ "/// " }), i(1, "Short description."),
    t({ "", "///", "/// " }), i(2, "Longer explanation (optional)."),
    t({ "", "///", "/// # Examples", "/// ```" }),
    t({ "", "/// # use rtc_rust::tuple::Tuple;" }),
    t({ "", "///" }),
    t({ "", "/// let result = " }), i(3, "Tuple::new(...)"),
    t({ ";", "///" }),
    t({ "", "/// assert!(" }), i(4, "true"), t({ ");" }),
    t({ "", "/// ```" }),
    t({ "", "" }),
  }),
})
