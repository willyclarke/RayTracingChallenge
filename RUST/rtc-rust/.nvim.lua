-- .nvim.lua — Rust project local config
-- This file is trusted per-project; load with :set exrc

vim.g.rustaceanvim = {
  server = {
    default_settings = {
      ['rust-analyzer'] = { cargo = { features = 'all' } },
    },
  },
}

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
local f = ls.function_node

-- <C-K> to expand (avoid <C-Y>, which conflicts with completion-menu accept)
vim.keymap.set({ "i" }, "<C-K>", function() ls.expand() end, { silent = true })
vim.keymap.set({ "i", "s" }, "<C-J>", function() ls.jump(-1) end, { silent = true })

vim.keymap.set({ "i", "s" }, "<C-E>", function()
  if ls.choice_active() then
    ls.change_choice(1)
  end
end, { silent = true })

-- <Tab>: jump through snippet placeholders, otherwise fall back to a real Tab
vim.keymap.set({ "i", "s" }, "<Tab>", function()
  if ls.expand_or_jumpable() then
    ls.expand_or_jump()
  else
    vim.api.nvim_feedkeys(
      vim.api.nvim_replace_termcodes("<Tab>", true, false, true),
      "n",
      false
    )
  end
end, { silent = true })

ls.add_snippets("rust", {
  s("tst", {
    t({ "/// Chap x - " }),
    i(1),                                  -- description: cursor lands here, paste from clipboard
    t({ "", "#[test]", "fn " }),
    i(2, "test_chap_x_y"),
    t({ "() -> Result<(), String> {",
        "    " }),
    i(3),
    t({ "",
        "    let chk = Tuple::point(1.0, 1.0, 1.0).approx_eq(Tuple::point(1.0, 1.0, 1.0));",
        "    if " }),
    i(4, "chk"),
    t({ " {",
        "        Ok(())",
        "    } else {",
        "        Err(\"" }),
    f(function(args) return args[1] end, { 1 }),  -- mirror i(1) into the Err string
    t({ "\".into())",
        "    }",
        "}" }),
    i(0),
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
